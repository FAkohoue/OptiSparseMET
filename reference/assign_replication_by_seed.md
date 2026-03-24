# Classify treatments into replication roles based on seed availability

Seed availability is the binding constraint that determines whether a
target replication level is deployable for any given treatment. This
function makes that constraint explicit by checking each treatment
against a per-plot seed requirement, an optional minimum buffer, and the
requested replication level, then assigning roles accordingly.
Treatments that cannot support the target replication may be downgraded
to unreplicated or excluded from the design, depending on
`shortage_action`.

Three design modes are supported. In `"augmented"` mode, all treatments
are unreplicated; the function identifies which ones have sufficient
seed for a single plot and flags the rest for exclusion. In `"p_rep"`
mode, a subset of treatments is selected for replication at
`desired_replications` and the remainder are unreplicated; selection is
governed by `priority`, restricted to `candidate_prep`, and optionally
capped at `max_prep`. In `"rcbd_type"` mode, all treatments are targeted
for `desired_replications`; those with insufficient seed are handled
according to `shortage_action`.

## Usage

``` r
assign_replication_by_seed(
  treatments,
  check_treatments = NULL,
  seed_available,
  seed_required_per_plot,
  replication_mode = c("augmented", "p_rep", "rcbd_type"),
  desired_replications = 2,
  candidate_prep = NULL,
  priority = c("seed_available", "input_order", "random"),
  minimum_seed_buffer = 0,
  max_prep = NULL,
  shortage_action = c("error", "downgrade", "exclude"),
  seed = NULL
)
```

## Arguments

- treatments:

  Character vector of non-check treatment IDs assigned to the
  environment. Check treatments should not appear here; they are managed
  separately within the field design functions. Duplicate values are
  silently deduplicated.

- check_treatments:

  Optional character vector of check IDs. If supplied, the function
  verifies that no treatment in `treatments` is also listed as a check,
  stopping with an error if overlap is found.

- seed_available:

  Data frame with columns `Treatment` and `SeedAvailable`.
  `SeedAvailable` is the total seed quantity available for each
  treatment. Every treatment in `treatments` must have a matching row;
  the function stops if any are missing.

- seed_required_per_plot:

  Positive numeric scalar. Seed quantity consumed by a single plot of
  one treatment in this environment. Used together with
  `desired_replications` and `minimum_seed_buffer` to compute
  feasibility thresholds.

- replication_mode:

  Character. Design mode controlling how the function partitions
  treatments into roles. One of:

  - `"augmented"`: targets one unreplicated plot per treatment.

  - `"p_rep"`: targets replication for a feasible subset, unreplicated
    for the remainder.

  - `"rcbd_type"`: targets `desired_replications` plots for all
    treatments.

  The selected mode must be consistent with the field design intended
  for the environment.
  [`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)
  interprets the output roles directly.

- desired_replications:

  Positive integer. Number of plots requested for replicated treatments.
  In `"augmented"` mode this value is not used for role assignment but
  is stored for reference. In `"p_rep"` mode it applies to the selected
  replicated subset. In `"rcbd_type"` mode it applies to all feasible
  treatments.

- candidate_prep:

  Optional character vector. Treatments eligible for replication in
  `"p_rep"` mode. Must be a subset of `treatments`; values not present
  in `treatments` are silently dropped. If `NULL`, all treatments are
  candidates. Ignored in `"augmented"` and `"rcbd_type"` modes.

- priority:

  Character. Selection rule for determining which candidate treatments
  are replicated in `"p_rep"` mode when the feasible candidate pool
  exceeds `max_prep`. One of:

  - `"seed_available"`: ranks by descending seed quantity, breaking ties
    by position in `treatments`.

  - `"input_order"`: follows the order of `treatments` as supplied.

  - `"random"`: draws uniformly at random from feasible candidates.

  Ignored in `"augmented"` and `"rcbd_type"` modes.

- minimum_seed_buffer:

  Non-negative numeric scalar. Additional seed quantity subtracted from
  the available supply before feasibility is assessed. Use this when a
  safety reserve must be maintained, for example to allow for
  germination losses or re-planting. The buffer is applied uniformly
  across all treatments.

- max_prep:

  Optional non-negative integer. Maximum number of treatments assigned
  the `"p_rep"` role in `"p_rep"` mode. If `NULL`, all treatments
  satisfying the replicated feasibility condition and present in
  `candidate_prep` are replicated. Use this when the field layout
  imposes an upper bound on the number of replicated plots independent
  of seed availability. Ignored in `"augmented"` and `"rcbd_type"`
  modes.

- shortage_action:

  Character. Action taken when a treatment has insufficient seed for its
  target role. One of:

  - `"error"`: stops immediately with an informative message listing the
    affected treatments.

  - `"downgrade"`: reduces the role to `"unreplicated"` if the
    single-plot feasibility condition is met; otherwise assigns
    `"excluded"`. Applies in `"rcbd_type"` mode only.

  - `"exclude"`: removes the treatment from the design.

- seed:

  Optional integer. Random seed for reproducibility. Affects random
  candidate selection in `"p_rep"` mode when `priority = "random"`. If
  `NULL`, no seed is set and results may differ across runs. The seed
  used internally is returned as `seed_used`.

## Value

A named list with the following components:

- `p_rep_treatments`:

  Character vector of treatments assigned the `"p_rep"` role, in the
  order they appear in `treatments`. Empty if no treatments are
  replicated.

- `p_rep_reps`:

  Integer vector of the same length as `p_rep_treatments` giving the
  replication count for each. All values equal `desired_replications`.

- `unreplicated_treatments`:

  Character vector of treatments assigned exactly one plot, in the order
  they appear in `treatments`.

- `excluded_treatments`:

  Character vector of treatments removed from the design due to
  insufficient seed, in the order they appear in `treatments`.

- `seed_summary`:

  Data frame with one row per treatment in `treatments`. Contains
  columns `Treatment`, `SeedAvailable`, `Required_Unrep` (seed threshold
  for one plot), `Required_Rep` (seed threshold for
  `desired_replications` plots), `Can_Unrep` (logical), `Can_Rep`
  (logical), and `Role` (the assigned role string).

- `replication_mode`:

  Character scalar. The mode used, returned for traceability.

- `seed_used`:

  The integer seed passed to the function, or `NULL` if none was
  supplied.

## Details

`assign_replication_by_seed()` evaluates seed feasibility for each
non-check treatment assigned to a single environment and partitions
those treatments into replication roles — replicated, unreplicated, or
excluded — according to the requested design mode and the available seed
quantity per treatment. The output is intended as direct input to
[`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md),
which constructs the within-environment field layout from the resulting
role assignments.

Let \\s_i\\ be the seed available for treatment \\i\\, \\q\\ the seed
required per plot, \\r\\ the requested replication count
(`desired_replications`), and \\b\\ the minimum seed buffer
(`minimum_seed_buffer`). The feasibility conditions are:

\$\$s_i \geq q + b \quad \text{(feasible for one unreplicated plot)}\$\$
\$\$s_i \geq r \times q + b \quad \text{(feasible for } r \text{
replicated plots)}\$\$

A treatment that fails the first condition cannot be placed in any plot
and is excluded regardless of mode. A treatment that satisfies the first
but not the second condition can support one plot but not the target
replication; its fate under shortage depends on `replication_mode` and
`shortage_action`.

### Replication modes

**`"augmented"`**: All treatments are targeted for a single unreplicated
plot. Treatments satisfying \\s_i \geq q + b\\ are assigned the role
`"unreplicated"`. Those that fail are assigned `"excluded"` if
`shortage_action` is `"exclude"` or `"downgrade"`, or trigger an error
if `shortage_action = "error"`.

**`"p_rep"`**: A selection of treatments is designated for replication
at `desired_replications`; the remaining treatments are unreplicated.
Candidates for replication are drawn from `candidate_prep` (all
treatments if `NULL`) restricted to those satisfying the replicated
feasibility condition. Among those candidates, selection priority is
determined by `priority`: `"seed_available"` ranks by descending seed
quantity, `"input_order"` preserves the order of `treatments`, and
`"random"` draws uniformly. The number of replicated treatments is
optionally capped by `max_prep`. Treatments not selected for replication
are evaluated for single-plot feasibility and assigned `"unreplicated"`
or `"excluded"`.

**`"rcbd_type"`**: All treatments are targeted for
`desired_replications`. Treatments satisfying the replicated feasibility
condition receive role `"p_rep"` at the requested replication level. For
treatments that fail, behavior depends on `shortage_action`: `"error"`
stops immediately; `"exclude"` drops the treatment from the design;
`"downgrade"` reduces the treatment to a single unreplicated plot if the
single-plot feasibility condition is met, or excludes it otherwise.

### Role vocabulary

The function assigns each treatment one of four roles, stored in the
`Role` column of `seed_summary`:

- `"p_rep"`: assigned `desired_replications` plots.

- `"unreplicated"`: assigned exactly one plot.

- `"excluded"`: assigned no plots; insufficient seed.

- `"unused"`: treatment not matched in the seed data frame (internal
  flag).

## Examples

``` r
treatments <- paste0("L", sprintf("%03d", 1:10))

seed_df <- data.frame(
  Treatment     = treatments,
  SeedAvailable = c(50, 45, 40, 35, 30, 20, 18, 15, 12, 10),
  stringsAsFactors = FALSE
)

## Example 1: p_rep mode — replicate the 4 best-seeded candidates
## at 2 reps; remaining feasible treatments receive one plot.
out1 <- assign_replication_by_seed(
  treatments             = treatments,
  seed_available         = seed_df,
  seed_required_per_plot = 10,
  replication_mode       = "p_rep",
  desired_replications   = 2,
  priority               = "seed_available",
  max_prep               = 4
)

out1$p_rep_treatments       # L001–L004: sufficient seed for 2 plots each
#> [1] "L001" "L002" "L003" "L004"
out1$unreplicated_treatments # remaining feasible treatments
#> [1] "L005" "L006" "L007" "L008" "L009" "L010"
out1$excluded_treatments     # treatments with seed < 10
#> character(0)
out1$seed_summary
#>    Treatment SeedAvailable Required_Unrep Required_Rep Can_Unrep Can_Rep
#> 1       L001            50             10           20      TRUE    TRUE
#> 2       L002            45             10           20      TRUE    TRUE
#> 3       L003            40             10           20      TRUE    TRUE
#> 4       L004            35             10           20      TRUE    TRUE
#> 5       L005            30             10           20      TRUE    TRUE
#> 6       L006            20             10           20      TRUE    TRUE
#> 7       L007            18             10           20      TRUE   FALSE
#> 8       L008            15             10           20      TRUE   FALSE
#> 9       L009            12             10           20      TRUE   FALSE
#> 10      L010            10             10           20      TRUE   FALSE
#>            Role
#> 1         p_rep
#> 2         p_rep
#> 3         p_rep
#> 4         p_rep
#> 5  unreplicated
#> 6  unreplicated
#> 7  unreplicated
#> 8  unreplicated
#> 9  unreplicated
#> 10 unreplicated

## Example 2: rcbd_type mode with downgrade — all treatments targeted for
## 3 reps; those below threshold are downgraded to 1 rep or excluded.
out2 <- assign_replication_by_seed(
  treatments             = treatments,
  seed_available         = seed_df,
  seed_required_per_plot = 10,
  replication_mode       = "rcbd_type",
  desired_replications   = 3,
  shortage_action        = "downgrade"
)

out2$p_rep_treatments        # treatments with seed >= 30
#> [1] "L001" "L002" "L003" "L004" "L005"
out2$unreplicated_treatments # treatments with 10 <= seed < 30
#> [1] "L006" "L007" "L008" "L009" "L010"
out2$excluded_treatments     # treatments with seed < 10
#> character(0)
out2$seed_summary
#>    Treatment SeedAvailable Required_Unrep Required_Rep Can_Unrep Can_Rep
#> 1       L001            50             10           30      TRUE    TRUE
#> 2       L002            45             10           30      TRUE    TRUE
#> 3       L003            40             10           30      TRUE    TRUE
#> 4       L004            35             10           30      TRUE    TRUE
#> 5       L005            30             10           30      TRUE    TRUE
#> 6       L006            20             10           30      TRUE   FALSE
#> 7       L007            18             10           30      TRUE   FALSE
#> 8       L008            15             10           30      TRUE   FALSE
#> 9       L009            12             10           30      TRUE   FALSE
#> 10      L010            10             10           30      TRUE   FALSE
#>            Role
#> 1         p_rep
#> 2         p_rep
#> 3         p_rep
#> 4         p_rep
#> 5         p_rep
#> 6  unreplicated
#> 7  unreplicated
#> 8  unreplicated
#> 9  unreplicated
#> 10 unreplicated

## Example 3: augmented mode with a seed buffer — each treatment gets one
## plot; a buffer of 5 seeds is reserved per treatment before feasibility
## is assessed.
out3 <- assign_replication_by_seed(
  treatments             = treatments,
  seed_available         = seed_df,
  seed_required_per_plot = 10,
  replication_mode       = "augmented",
  minimum_seed_buffer    = 5
)

out3$unreplicated_treatments # treatments with seed >= 15
#> [1] "L001" "L002" "L003" "L004" "L005" "L006" "L007" "L008"
out3$excluded_treatments     # treatments with seed < 15
#> [1] "L009" "L010"
out3$seed_summary
#>    Treatment SeedAvailable Required_Unrep Required_Rep Can_Unrep Can_Rep
#> 1       L001            50             15           25      TRUE    TRUE
#> 2       L002            45             15           25      TRUE    TRUE
#> 3       L003            40             15           25      TRUE    TRUE
#> 4       L004            35             15           25      TRUE    TRUE
#> 5       L005            30             15           25      TRUE    TRUE
#> 6       L006            20             15           25      TRUE   FALSE
#> 7       L007            18             15           25      TRUE   FALSE
#> 8       L008            15             15           25      TRUE   FALSE
#> 9       L009            12             15           25     FALSE   FALSE
#> 10      L010            10             15           25     FALSE   FALSE
#>            Role
#> 1  unreplicated
#> 2  unreplicated
#> 3  unreplicated
#> 4  unreplicated
#> 5  unreplicated
#> 6  unreplicated
#> 7  unreplicated
#> 8  unreplicated
#> 9      excluded
#> 10     excluded
```
