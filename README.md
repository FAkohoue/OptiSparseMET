# OptiSparseMET


[![R-CMD-check](https://github.com/FAkohoue/OptiSparseMET/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FAkohoue/OptiSparseMET/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

> Sparse multi-environment trial (MET) design that couples **across-environment
> allocation** with **within-environment field layout** under shared genetic,
> environmental, and seed constraints.

---

## Overview

`OptiSparseMET` is an R framework for designing sparse METs. Modern breeding
programmes evaluate more candidate genotypes than any one environment can hold,
environments differ in capacity and precision, and finite seed stocks constrain
both presence and replication. Deciding these separately can yield a design that
is statistically attractive but unplantable, or a feasible field book with weak
cross-environment information. `OptiSparseMET` coordinates all of it —
environmental evidence, common + sparse allocation, network-wide seed
accounting, local field layout, and robust pre-deployment evaluation — and
optimises the two design levels *together*.

The package distinguishes evidence from assumptions: historical MET responses,
when available, calibrate the genetic environmental covariance; when they are
not, it uses a neutral central covariance plus explicit sensitivity scenarios
rather than asking you to invent weights.

At its core it optimises the reliability of across-environment breeding values,
whose precision comes from a **two-level information matrix** that couples the
allocation (through the per-cell information $D$) with the environmental
covariance $\Sigma_E$ and the genomic relationships $G$:

$$
C_{uu} \;=\; D \;-\; D X (X^\top D X)^{-1} X^\top D \;+\;
\sigma_g^{-2}\,\big(\Sigma_E^{-1} \otimes G^{-1}\big).
$$

📘 **[Download the Breeder's Guide (PDF)](https://github.com/FAkohoue/OptiSparseMET/raw/master/inst/guides/OptiSparseMET_Breeders_Guide.pdf)** ·
📖 **[Full documentation & tutorials](https://FAkohoue.github.io/OptiSparseMET/)**

<p align="center">
  <img src="man/figures/OptiSparseMET_schematic.png"
       alt="OptiSparseMET workflow schematic" width="100%">
</p>

---

## What it does

The pipeline runs in nine modules (86 exported functions), each documented in
the [reference index](https://FAkohoue.github.io/OptiSparseMET/reference/) and
demonstrated end-to-end in the pipeline vignette:

1. **Environmental classification** — weather/soil/management kernels, a
   data-driven environmental covariance `Sigma_E`, and mega-environments.
2. **Genetic relationship matrices** — genomic, hybrid/testcross (with optional
   dominance), and made-invertible relationship matrices.
3. **Sparse allocation** — M3 (random-balanced) / M4 (equireplicate) allocation
   with feasibility checks and auto-suggested common treatments.
4. **Seed-aware replication** — a single seed inventory turned into a feasible
   per-site replication plan.
5. **Within-environment field design** — block and alpha row-column layouts with
   efficiency evaluation.
6. **Coupled optimisation & information** — the two-level information matrix,
   exact multi-trait index reliability, robust/CVaR optimisation, and Pareto
   frontiers.
7. **Simulation** — Monte-Carlo prediction accuracy and realised genetic gain.
8. **Benchmarking** — paired comparison of reference designs with confidence
   intervals, tail risk, and decision stability.
9. **Field books** — per-site and combined MET field books, and field maps.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("FAkohoue/OptiSparseMET",
                        build_vignettes = TRUE, dependencies = TRUE)
```

Set `build_vignettes = FALSE` for a faster install.

---

## Quick start

```r
library(OptiSparseMET)

# Operational one-call pipeline: allocation -> seed -> local design -> field book
out <- plan_sparse_met_design(
  treatments                     = sprintf("H%03d", 1:120),
  environments                   = c("E1", "E2", "E3", "E4"),
  allocation_method              = "random_balanced",
  n_test_entries_per_environment = 41,
  target_replications            = 1
)
out$combined_field_book   # the assembled MET field book
```

For the **full scientific pipeline** — from environmental classification through
optimisation, simulation, and benchmarking to the combined field book — see the
tutorial:

```r
vignette("OptiSparseMET-pipeline", package = "OptiSparseMET")
```

---

## Documentation

| Resource | Contents |
|----------|----------|
| [**Pipeline tutorial**](https://FAkohoue.github.io/OptiSparseMET/articles/OptiSparseMET-pipeline.html) | Full end-to-end run of all nine modules on a reproducible example |
| [Introduction](https://FAkohoue.github.io/OptiSparseMET/articles/OptiSparseMET-introduction.html) | Statistical framework, feasibility rules, and input contract |
| [Environmental interactions](https://FAkohoue.github.io/OptiSparseMET/articles/OptiSparseMET-environmental-interactions.html) | Enviromic kernels and interaction evidence |
| [Benchmarking](https://FAkohoue.github.io/OptiSparseMET/articles/OptiSparseMET-benchmarking.html) | Comparing and validating designs before release |
| [Breeder's Guide (PDF)](https://github.com/FAkohoue/OptiSparseMET/raw/master/inst/guides/OptiSparseMET_Breeders_Guide.pdf) | Decision-oriented companion, no R required |
| [Function reference](https://FAkohoue.github.io/OptiSparseMET/reference/) | All 86 functions, grouped by module |

```r
# After installation:
vignette(package = "OptiSparseMET")   # list all vignettes
```

---

## Citation

If you use `OptiSparseMET` in published research, please cite:

```
Akohoue, F. (2026).
OptiSparseMET: Sparse Multi-Environment Trial Design with Flexible Local
Field Layout. R package version 0.2.0.
https://github.com/FAkohoue/OptiSparseMET
```

## Reference

Montesinos-Lopez O.A., Mosqueda-Gonzalez B.A., Salinas-Ruiz J.,
Montesinos-Lopez A., Crossa J. (2023). Sparse multi-trait genomic prediction
under balanced incomplete block design. *The Plant Genome*, 16, e20305.
<https://doi.org/10.1002/tpg2.20305>

## Contributing

Issues, bug reports, and feature suggestions are welcome:
<https://github.com/FAkohoue/OptiSparseMET/issues>

## License

MIT License © Félicien Akohoue
