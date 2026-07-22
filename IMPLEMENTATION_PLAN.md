# OptiSparseMET — Implementation plan to expert-grade

Prepared for F. Akohoue — 2026-07-21. Scope: **Tiers 1–3** (correctness + rigor,
validation + coupling, and the Colmant-framework expansion). Companion to
`CRITIQUE_Colmant2026_and_OptiSparseMET.md` — each item below cites the critique
tag it closes (e.g. B1a, C3).

Each work item has: **Problem → Files → Change → Test (with oracle) → Done**.
"Oracle" = the independent source of truth the test checks against (closed-form
value, hand-built matrix, or reference package), never the code's own formula.

Legend for dependencies: items are grouped into phases P0–P6 in execution order.
Later phases assume earlier ones are merged.

---

## Phase P0 — Correctness bugs (do first; highest trust impact)

### P0.1 — AR1 / AR1×AR1 residual Kronecker indexing (closes B1a)
- **Problem:** `grid_index = (Row-1)*n_cols + Column` indexes `kronecker(Qcol,
  Qrow)` as if it were `kronecker(Qrow, Qcol)`. Row/column autocorrelation axes
  are swapped and adjacency scrambles on non-square fields.
- **Files:** `R/met_evaluate_alpha_efficiency.R` (~L376–388),
  `R/met_evaluate_famoptg_efficiency.R` (~L380–391). Same bug in both.
- **Change:** either (a) keep `kronecker(Qcol, Qrow)` and set
  `grid_index = (Column-1)*n_rows + Row`, or (b) keep the current index and
  build `kronecker(Qrow, Qcol)`. Pick (a) for both files; add an inline comment
  deriving the index. Do the same fix for the row-only `AR1` branch
  (`kronecker(Diagonal(n_cols), Qrow)` → index `(Column-1)*n_rows + Row`).
- **Test** (`tests/testthat/test-ar1-precision.R`, new): build a **non-square**
  3×4 field with `rho_row = 0.5`, `rho_col = 0.2`. Oracle: form
  `Sigma_row[i,j] = 0.5^|i-j|` (3×3), `Sigma_col = 0.2^|i-j|` (4×4),
  `Sigma = kronecker(Sigma_col, Sigma_row)`, `Q_oracle = solve(Sigma)`. Map each
  field plot to its oracle index and assert the code's `Q` (×`sigma_e2`) equals
  `Q_oracle` up to `1e-10`. Add a guard assertion that swapping `rho_row`↔`rho_col`
  changes `A_criterion` (proves axes are distinct). Regression: this test fails
  on current code, passes after the fix.
- **Done:** both evaluators reproduce the hand-built AR1×AR1 precision on a
  non-square grid; SA/GA runs under `residual_structure="AR1xAR1"` unaffected in
  API.

### P0.2 — CDmean ignores diag(K) (closes B1b / C2)
- **Problem:** `CDmean = 1 - mean_PEV/sigma_g2` assumes `K_ii = 1`; wrong for
  VanRaden GRM / inbreds (diag ≈ 1+F).
- **Files:** `R/met_evaluate_alpha_efficiency.R` (random branch, ~L601–648),
  `R/met_evaluate_famoptg_efficiency.R` (matching block).
- **Change:** compute per-line `CD_i = 1 - PEV_i / (sigma_g2 * Kdiag_i)` where
  `Kdiag_i = diag(K)` aligned to the line order used for `Zg` (reuse the existing
  `ids`/`line_ids` mapping; for `prediction_type="IID"` set `Kdiag_i = 1`).
  `CDmean = mean(CD_i)`. Keep `mean_PEV` reported as-is. Add a `cd_definition`
  field ("per_line_reliability") and update the roxygen to stop citing Rincent
  as if it were the contrast-based CDmean, or implement the true contrast form
  as `CDmean_contrast` (optional, see P3.4).
- **Test** (`test-cdmean.R`, new): (1) backward-compat — with `diag(K)=1`, new
  `CDmean` equals old `1 - mean_PEV/sigma_g2` to `1e-10`. (2) Correctness — set
  `diag(K)=1.5`; oracle computes `PEV_i` from an independent small MME solve and
  asserts `CD_i == 1 - PEV_i/(sigma_g2*1.5)` and all `CD_i ∈ [0,1]`. (3) Guard —
  a case where the old formula produced `CDmean > 1` now stays in range.
- **Done:** CDmean correct for non-unit-diagonal K; docs no longer mislabel it.

### P0.3 — D-criterion pseudo-determinant tolerance (closes C2)
- **Problem:** `.safe_logdet_psd_dense` uses absolute `tol = 1e-10`, but the
  divisor `q = p-1` is hardcoded; eigenvalue count and divisor can disagree.
- **Files:** `R/alpha_rc_helpers.R` (`.safe_logdet_psd_dense`); callers in both
  evaluators (pass expected rank).
- **Change:** relative tolerance `tol_rel * max(ev)` (default `tol_rel=1e-8`);
  optionally accept an `expected_rank` arg and, after thresholding, keep exactly
  the top `expected_rank` eigenvalues (here `p-1`), erroring if the count is off
  by more than rounding. Return both `logdet` and `n_pos` so the caller can
  assert `n_pos == q`.
- **Test** (`test-logdet.R`, new): build `V = diag(c(4, 1, 1e-9))` scale, centered
  `HVH`; oracle geometric mean over the intended `p-1` eigenvalues. Assert
  `D_criterion` matches oracle and is invariant to multiplying V by a constant ×
  then dividing (scale robustness). Assert a near-zero structural eigenvalue is
  dropped and a genuine small eigenvalue is retained.
- **Done:** D is scale-robust and always averages exactly `p-1` eigenvalues.

### P0.4 — Rank-deficiency guard in fixed branch (closes C2)
- **Problem:** `check_as_fixed=FALSE` + `treatment_effect="fixed"` keeps the
  intercept alongside a full dummy set → singular C → `.solve_C` errors.
- **Files:** `R/met_evaluate_alpha_efficiency.R`,
  `R/met_evaluate_famoptg_efficiency.R` (fixed-effects assembly).
- **Change:** when the fixed design is rank-deficient, drop the intercept (or a
  reference level) and emit an informative message; detect via
  `Matrix::rankMatrix` or a QR pivot on `XtQX`.
- **Test** (`test-rank-guard.R`, new): call with `check_as_fixed=FALSE`,
  `treatment_effect="fixed"`; assert it returns finite `A_criterion` (no error)
  and that the value equals the `check_as_fixed=TRUE` result on an equivalent
  parameterization (contrast-invariance).
- **Done:** no singular-matrix crash; contrasts unaffected.

---

## Phase P1 — Allocation: relabel M4 and deliver achievable balance (closes C3)

### P1.1 — Rename the BIBD claim honestly
- **Problem:** `balanced_incomplete`/`check_balanced_incomplete_feasibility()`
  promise a BIBD (constant λ) that is unachievable when J ≫ I and is not
  constructed.
- **Files:** `R/allocate_sparse_met.R`, `R/check_balanced_incomplete_feasibility.R`,
  `NAMESPACE`, `man/*`, `README.md`, `DESCRIPTION`.
- **Change:** keep `"balanced_incomplete"` as a **deprecated alias** but
  introduce `allocation_method = "equireplicate"` (equal r + equal k) as the
  documented name. Add `.Deprecated`-style messaging. Reword docs: "equireplicate
  incomplete-block allocation; approaches a near-balanced (regular-graph) design;
  a strict BIBD generally does not exist for J ≫ I (Fisher's inequality)."
- **Test** (`test-allocate_naming.R`): assert the old name still works and emits
  a deprecation warning; new name produces identical structure for a seed-fixed
  case.
- **Done:** naming matches what is delivered; back-compat preserved.

### P1.2 — Near-balanced (regular-graph) line-pair concurrence
- **Problem:** the greedy constructor equalizes r and k but never balances
  pairwise concurrence, the property that drives connectivity.
- **Files:** `R/allocate_sparse_met.R` (strict constructor 7A), new helper
  `R/allocation_balance_helpers.R`.
- **Change:** after the equireplicate assignment, run a **swap-improvement pass**:
  target average concurrence `lambda_bar = r*(k_star-1)/(J_star-1)`; minimize
  `sum((overlap_offdiag - lambda_bar)^2)` by proposing treatment↔treatment swaps
  between environments that preserve r and k, accepting swaps that reduce the
  objective (optional SA wrapper). Expose `balance_concurrence = TRUE` and
  `balance_iter`.
- **Test** (`test-concurrence-balance.R`): for parameters where a regular-graph
  design exists (e.g. J*=9, I=3, k*=3, r=1 → all λ∈{0,1}; and a case with
  λ_bar≈1.x), assert off-diagonal concurrences take only the two adjacent values
  `{floor(lambda_bar), ceil(lambda_bar)}` and that variance of concurrences is
  ≤ the pre-balance variance. Oracle: exact two-value set computed from
  `lambda_bar`.
- **Done:** allocation is a documented near-balanced design; concurrence variance
  minimized.

### P1.3 — Environment-pair balance (cross-environment connectivity)
- **Problem:** the more meaningful MET property (equal shared-line counts between
  environment pairs) is not targeted.
- **Files:** `R/allocation_balance_helpers.R`, `R/allocate_sparse_met.R`.
- **Change:** compute the environment-by-environment shared-line matrix
  `E_share = alloc_sparse^T alloc_sparse` off-diagonal; add option
  `balance = c("line_pair","env_pair","both")` with a weighted objective when
  `"both"` (document the trade-off — they can conflict). For `"env_pair"`,
  minimize variance of `E_share` off-diagonals via the same swap engine.
- **Test** (`test-env-pair-balance.R`): for J*=100, I=5, r=2, assert env-pair
  shared-line counts fall within ±1 of their mean after balancing, and that
  `balance="env_pair"` reduces `var(E_share_offdiag)` vs unbalanced. Oracle:
  expected mean shared count `= r*(r-1)/(I-1) * ...` computed directly from the
  incidence totals.
- **Done:** users can prioritize the connectivity property MET inference depends
  on.

### P1.4 — Real feasibility conditions
- **Problem:** `check_balanced_incomplete_feasibility()` only checks the slot
  identity, not λ-integrality or Fisher's inequality.
- **Files:** `R/check_balanced_incomplete_feasibility.R`.
- **Change:** rename to `check_equireplicate_feasibility()` (keep old as alias);
  return `slot_ok`, plus informational `bibd_lambda = r*(k_star-1)/(J_star-1)`,
  `bibd_lambda_integer` (logical), `fisher_ok = I >= J_star`, and a
  `strict_bibd_possible` flag = `slot_ok & lambda_integer & fisher_ok`. Do **not**
  block on the BIBD conditions (they are informational).
- **Test** (`test-feasibility.R`, extend): assert `strict_bibd_possible=FALSE`
  for J*=110,I=4,r=2 (Fisher fails) and `TRUE` for a Fano-type v=b=7,r=k=3,λ=1.
- **Done:** the checker reports what is actually true about BIBD existence.

### P1.5 — Optional exact BIBD via crossdes when it exists
- **Problem:** in the small-J regime a true BIBD exists but the greedy
  constructor won't find it; the paper uses `crossdes::find.BIB()`.
- **Files:** `R/allocate_sparse_met.R`; `DESCRIPTION` (Suggests: crossdes).
- **Change:** when `strict_bibd_possible` and `crossdes` is installed, offer
  `construction = "exact_bibd"` delegating to `crossdes::find.BIB()` +
  `crossdes::isBIB()` verification; else fall back with a message.
- **Test** (`test-exact-bibd.R`, `skip_if_not_installed("crossdes")`): construct
  the v=b=7 BIBD, assert `isBIB()` TRUE and constant λ=1.
- **Done:** exact BIBD available where mathematically possible.

---

## Phase P2 — Replace tautological tests with ground-truth tests (closes B1c)

### P2.1 — Closed-form A/D oracles
- **Files:** `tests/testthat/test-efficiency-groundtruth.R` (new).
- **Change/Test:** CRD with `t` treatments each replicated `n`, IID residual
  `sigma_e2`, no blocking. Oracle: `Var(tau_i - tau_j) = 2*sigma_e2/n` ⇒
  `A_criterion == 2*sigma_e2/n`; cell-mean covariance `= sigma_e2/n * I` ⇒
  `D_criterion == sigma_e2/n`. Assert both to `1e-8`. Add an RCBD case with a
  known blocking-efficiency factor.
- **Done:** A and D verified against analytic truth, not against themselves.

### P2.2 — PEV / CDmean cross-check vs sommer
- **Files:** `tests/testthat/test-pev-vs-sommer.R` (new,
  `skip_if_not_installed("sommer")`).
- **Change/Test:** simulate a small trial + GRM; fit GBLUP in `sommer`; extract
  its PEV; assert package `mean_PEV` matches within Monte-Carlo/numeric tolerance
  (e.g. 2%). Confirms the MME assembly and C⁻¹ extraction.
- **Done:** PEV path validated against an independent mixed-model engine.

### P2.3 — Purge self-referential assertions
- **Files:** `test-met_evaluate_alpha_efficiency.R`,
  `test-met_evaluate_famoptg_efficiency.R`.
- **Change:** remove/repurpose `expect_equal(A_efficiency, 1/A_criterion)` and
  `expect_equal(CDmean, 1 - mean_PEV/sigma_g2)`; keep only as trivial identity
  smoke-checks clearly labeled as such, not as correctness evidence.
- **Done:** the suite would now catch P0.1–P0.3 regressions.

---

## Phase P3 — Validation, coupling, and honest criteria (Tier 2)

### P3.1 — A/D "efficiency" naming (closes C1)
- **Files:** both evaluators; docs.
- **Change:** deprecate bare `$A`/`$D` aliases (keep with a one-time warning);
  document `A_criterion` (lower better) vs `A_inv = 1/A_criterion`; add a true
  **relative** `A_efficiency ∈ (0,1]` computed against an orthogonal/ideal
  reference design when one is well-defined, else return `NA` with a note.
- **Test:** assert relative `A_efficiency ∈ (0,1]` for a sub-optimal design and
  `== 1` (within tol) for the reference design.
- **Done:** "efficiency" means the textbook relative quantity.

### P3.2 — MET-level (across-environment) information matrix (closes B2g)
- **Problem:** the headline "coupled two-level" claim is not implemented;
  evaluation is per-environment, single-trait.
- **Files:** new `R/met_information.R`; wire into `plan_sparse_met_design()`.
- **Change:** assemble the joint MME across environments using the allocation
  incidence `Z` (which lines in which environments) and a G×E covariance
  `G ⊗ Sigma_E` (start diagonal `Sigma_E`, then unstructured/FA option), with the
  within-environment `R^{-1}` blocks from the local designs. Return an
  across-TPE `mean_PEV`/`CDmean` for the *combined* design — the quantity the
  allocation exists to optimize.
- **Test** (`test-met-information.R`): two-environment toy with known G, Sigma_E;
  oracle builds the full `C` densely and inverts; assert the module's PEV matches.
  Also assert that increasing cross-environment connectivity (more shared lines)
  lowers across-TPE mean PEV.
- **Done:** the two-level coupling is computed, not just asserted in prose.

### P3.3 — G×E simulation harness → accuracy & genetic gain (closes B2d)
- **Problem:** proxies never linked to realized accuracy; the vignette
  "simulation" has no G×E.
- **Files:** new `R/simulate_met.R` (or a `data-raw` + vignette module);
  rewrite the vignette simulation section.
- **Change:** simulate genotype values with a real across-environment genetic
  **correlation structure** (e.g. `g ~ MVN(0, G ⊗ Sigma_E)`), add spatial +
  residual noise, fit GBLUP, and report accuracy (cor(true, predicted BV)) and
  Δg (top-10% selection) for designs the package produces. Provide a
  `benchmark_designs()` wrapper comparing balanced vs sparse vs near-balanced.
- **Test** (`test-simulate-met.R`, light): deterministic-seed smoke test that
  accuracy is finite, in [-1,1], and that a fully balanced design beats a
  degenerate family-clustered design (sanity direction, matching Colmant Fig 4).
- **Done:** the package can demonstrate design → accuracy/gain, not just A/D/CD.

### P3.4 — True Rincent contrast-based CDmean (optional, complements P0.2)
- **Files:** evaluators.
- **Change:** add `CDmean_contrast` using `M = I - X(X'X)^{-1}X'` and the
  generalized `CD(c)` normalized by `c'Kc`, so the Rincent citation is earned.
- **Test:** on a design with a fixed environment effect, assert
  `CDmean_contrast` is invariant to adding a constant to the fixed mean (contrast
  property) whereas raw reliability is not.
- **Done:** both the simple reliability and the true Rincent CDmean are available
  and correctly labeled.

### P3.5 — Variance-component sensitivity + optimizer diagnostics (closes B2e/B2f)
- **Files:** `R/met_optimize_alpha_rc.R`, `R/met_optimize_famoptg.R`, new
  `R/sensitivity.R`.
- **Change:** (a) `evaluate_over_varcomp_grid()` sweeping σ ratios / ρ and
  reporting whether the chosen design stays best (robustness table). (b) Return
  optimizer traces (score vs iteration/restart), restart-to-restart variance, and
  a `converged` flag; add a `benchmark = TRUE` mode comparing to a brute-force
  optimum on small problems.
- **Test:** assert the SA/GA best score ≤ random-restart best on a fixed small
  problem, and that the brute-force optimum is matched on a tiny enumerable case.
- **Done:** optimizers are diagnosable and shown near-optimal on small cases.

### P3.6 — External benchmarks (closes B3l)
- **Files:** `vignettes/benchmarks.Rmd` (new), `tests` behind `skip_if_not_installed`.
- **Change:** reproducible comparisons of A/D and PEV vs `odw`/`DiGGer`/`blocksdesign`
  on shared example inputs; short validation report.
- **Done:** published evidence the criteria/optimizers agree with reference tools.

---

## Phase P4 — Colmant decision-point framework (Tier 3)

### P4.1 — Environment selection from a TPE (−q′Dq + baselines)
- **Files:** new `R/select_environments.R`; export in `NAMESPACE`.
- **Change:** implement `select_environments(D, n, method = c("optcontrib",
  "kmeans","hclust","random"))` where `optcontrib` maximizes representativeness
  via `-q'Dq` under a cardinality constraint (GA or exact for small n), with the
  documented caveat that `-q'Dq` maximizes spread, not centroid-representativeness
  (offer an optional mean-representativeness term). Accept an environment
  relationship matrix `E` (from historical genetic correlations or enviromic data).
- **Test** (`test-select-env.R`): on a constructed `D` with a known redundant
  cluster, assert `optcontrib` avoids selecting two near-identical environments
  that `random` would pick; assert all methods honor the cardinality `n`.
- **Done:** decision points 1–2 available with honest documentation of the
  objective's meaning.

### P4.2 — TPG sparsity dimension
- **Files:** new `R/sparsity_grid.R`.
- **Change:** helper to explore the `ia × ipf × noe` trade-off under a total-plot
  constraint (Colmant §2.3), returning accuracy/Δg from the P3.3 harness for each
  cell so programs can locate their own sweet spot.
- **Test:** grid enumeration respects the plot budget; monotonicity sanity
  (accuracy non-decreasing in noe at fixed budget in the broad-TPE simulation).
- **Done:** the largest-effect decision in the paper is explorable in-package.

### P4.4 — Resource- and noise-aware replication recommendation (decision 4)
- **Files:** new `R/recommend_replication.R`, built on
  `R/assign_replication_by_seed.R`; uses the P3.3 harness.
- **Change:** given the program's spatial-variance level, GRM, and seed
  availability, sweep replication level × p-rep fraction × field size and return
  the accuracy/Δg-per-cost curve; recommend the diminishing-returns level.
  De-confound field size from unique-entry count (hold one fixed while varying
  the other), correcting Colmant §2.6. Never exceed seed limits.
- **Test** (`test-recommend-replication.R`): recommended replication differs
  between a low-noise and a high-noise scenario (noise-dependence the
  single-setting paper could not show); recommendations respect seed availability.
- **Done:** decision 4 answered by a resource-aware criterion, not a static grid.

### P4.3 — Objective-function-driven G×E allocation
- **Files:** `R/allocate_sparse_met.R` (new `objective = "optcontrib_GxE"`).
- **Change:** allow allocation to optimize `-q'Dq` with `D = G ⊗ E` (Colmant
  §2.5) via the swap/GA engine from P1, as an alternative to combinatorial
  equireplicate — unifying the paper's approach with the package's constraints
  (seed, capacity, common treatments).
- **Test:** assert the objective decreases monotonically over iterations and that
  the result improves across-TPE PEV (P3.2) vs random allocation on a broad-TPE
  case.
- **Done:** the paper's flagship allocation objective is available and validated
  against the coupling metric.

---

## Phase P5 — Cross-cutting: packaging, CI, docs

- **DESCRIPTION/NEWS:** bump to 0.2.0; document deprecations; add `crossdes`,
  `sommer` (Suggests). Update `Config/lifecycle` as components stabilize.
- **NAMESPACE:** export new functions (`select_environments`,
  `check_equireplicate_feasibility`, `simulate_met`, `benchmark_designs`, …).
- **CI:** ensure `R-CMD-check` runs the new suites; add `covr` coverage gate
  (target ≥ 85% on evaluators and allocator); add a `testthat` snapshot for
  design summaries.
- **README:** replace BIBD language with "equireplicate / near-balanced"; state
  clearly that strict BIBD is generally infeasible for J ≫ I; move the
  "coupled two-level" claim next to the P3.2 function that now backs it.
- **Reproducibility:** meaningful commit messages; tag `v0.2.0`; add a
  `VALIDATION.md` summarizing P2/P3.6 results.

---

## Suggested execution order (dependency-aware)

1. **P0.1–P0.4** (bugs) — independent, parallelizable, merge first.
2. **P2.1–P2.3** (ground-truth tests) — lock the corrected behavior.
3. **P1.1–P1.5** (allocation relabel + balance + feasibility).
4. **P3.1** (naming), **P3.2** (MET information matrix) — unblocks validation.
5. **P3.3** (simulation harness) — depends on P3.2 for the coupling metric.
6. **P3.4–P3.6** (Rincent CDmean, sensitivity/diagnostics, benchmarks).
7. **P4.1–P4.3** (framework expansion) — depends on P3.2/P3.3 harness.
8. **P5** — finalize packaging, CI, docs, validation report.

**Definition of done for the whole effort:** every criterion verified against an
external oracle (not itself); the allocator delivers and correctly names a
near-balanced design with reported concurrence/connectivity balance; the
across-TPE information matrix and a G×E simulation demonstrate that the
package's designs raise realized accuracy/genetic gain; and the framework
(environment selection, TPG sparsity, G×E allocation) is available with honest
documentation of each objective's assumptions.
