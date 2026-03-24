# OptiSparseMET: Sparse Multi-Environment Trial Design with Flexible Local Field Layout Construction

`OptiSparseMET` constructs sparse multi-environment trial designs by
addressing treatment allocation across environments and field layout
construction within environments as a two-level problem, solved jointly
under shared statistical and genetic constraints. The package is
intended for breeding programs where the number of candidate lines
exceeds what any single environment can accommodate, environments are
heterogeneous in capacity and genetic composition, and seed availability
imposes practical limits on replication that purely theoretical designs
ignore.

Package providing tools for:

- sparse treatment allocation across environments (M3 and M4
  strategies),

- feasibility checking and capacity helpers,

- seed-aware replication planning,

- within-environment field layout construction,

- MET-level field book assembly,

- end-to-end pipeline execution.

## Details

\*\* Stage 1: across-environment sparse allocation \*\*

[`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)
determines which treatments enter which environments. Two allocation
strategies are available. `"random_balanced"` implements an M3-type
stochastic allocation that approximates balance without requiring exact
BIBD parameters to be satisfiable — appropriate when environment
capacities differ or when the trial dimensions do not admit an exact
balanced solution. `"balanced_incomplete"` implements an M4-type
allocation following BIBD principles at the MET level, enforcing equal
replication and approximately uniform pairwise co-occurrence —
appropriate when environments are comparable in size and equal
replication is a hard requirement.

Allocation can be guided by genetic structure through family labels,
genomic relationship clusters (GRM), or pedigree relationship clusters
(A), ensuring that genetic groups are distributed across environments
rather than concentrated within a subset of them. Common treatments can
be forced into all environments to establish direct cross-environment
connectivity independent of model assumptions. Feasibility of an exact
balanced incomplete allocation can be checked in advance with
[`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md).

The underlying resource identity is:

\$\$N = J \times r = I \times k\$\$

where \\J\\ is the number of treatments, \\r\\ is the number of
environments each treatment enters, \\I\\ is the number of environments,
and \\k\\ is the number of treatments per environment. This identity
makes the tradeoff between coverage and replication depth explicit given
fixed total resources \\N\\.

\*\* Stage 2: within-environment field design \*\*

Each environment is designed independently after allocation.
[`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)
constructs augmented, partially replicated (p-rep), or RCBD-type
repeated-check designs for environments where block-based local control
is appropriate.
[`alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/alpha_rc_stream.md)
constructs stream-based alpha row-column designs for environments with
fixed grid geometry and continuous field-book order, accommodating AR1
and AR1×AR1 spatial models at the analysis stage.

When replication within an environment is flexible,
[`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)
evaluates seed availability against the per-plot seed requirement for
each treatment and assigns replication roles — replicated, unreplicated,
or excluded — before the field layout is constructed. Both
within-environment design functions accept optional relationship
matrices for structure-aware entry arrangement and dispersion
optimization, and can evaluate the resulting design under a mixed-model
framework prior to field deployment.

Environment-level field books are assembled into a single MET field book
with
[`combine_met_fieldbooks()`](https://FAkohoue.github.io/OptiSparseMET/reference/combine_met_fieldbooks.md),
which handles heterogeneous column sets across environments and adds
MET-level metadata to every row. The complete two-stage pipeline can
also be run end-to-end through
[`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md).

## Exported functions

|                                                                                                                                          |                                                       |
|------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------|
| Function                                                                                                                                 | Role                                                  |
| [`allocate_sparse_met()`](https://FAkohoue.github.io/OptiSparseMET/reference/allocate_sparse_met.md)                                     | Allocate treatments across environments               |
| [`check_balanced_incomplete_feasibility()`](https://FAkohoue.github.io/OptiSparseMET/reference/check_balanced_incomplete_feasibility.md) | Verify slot identity before allocation                |
| [`derive_allocation_groups()`](https://FAkohoue.github.io/OptiSparseMET/reference/derive_allocation_groups.md)                           | Derive genetic group labels for allocation            |
| [`assign_replication_by_seed()`](https://FAkohoue.github.io/OptiSparseMET/reference/assign_replication_by_seed.md)                       | Classify treatments into replication roles            |
| [`prep_famoptg()`](https://FAkohoue.github.io/OptiSparseMET/reference/prep_famoptg.md)                                                   | Build block-based within-environment layouts          |
| [`alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/alpha_rc_stream.md)                                             | Build stream-based row-column layouts                 |
| [`combine_met_fieldbooks()`](https://FAkohoue.github.io/OptiSparseMET/reference/combine_met_fieldbooks.md)                               | Stack environment field books into one MET field book |
| [`plan_sparse_met_design()`](https://FAkohoue.github.io/OptiSparseMET/reference/plan_sparse_met_design.md)                               | Run the full two-stage pipeline                       |

## Example data

`OptiSparseMET_example_data` is a synthetic dataset containing treatment
vectors, environment names, family metadata, seed inventories,
relationship matrices (GRM, A, K), and pre-built argument lists for
allocation and local design functions. All components refer to the same
set of synthetic lines and environments and can be passed to any package
function without modification. Load with:

    data("OptiSparseMET_example_data", package = "OptiSparseMET")

## Dependencies

**Matrix** is used for sparse matrix operations in efficiency evaluation
and design construction. **pracma** provides `mod()` for serpentine
traversal logic in
[`alpha_rc_stream()`](https://FAkohoue.github.io/OptiSparseMET/reference/alpha_rc_stream.md).

## References

Montesinos-López, O. A., Mosqueda-González, B. A., Salinas-Ruiz, J.,
Montesinos-López, A., & Crossa, J. (2023). Sparse multi-trait genomic
prediction under balanced incomplete block design. *The Plant Genome*,
16, e20305.

## See also

Useful links:

- <https://github.com/FAkohoue/OptiSparseMET>

- <https://FAkohoue.github.io/OptiSparseMET>

- Report bugs at <https://github.com/FAkohoue/OptiSparseMET/issues>

Useful links:

- <https://github.com/FAkohoue/OptiSparseMET>

- <https://FAkohoue.github.io/OptiSparseMET>

- Report bugs at <https://github.com/FAkohoue/OptiSparseMET/issues>

## Author

**Maintainer**: Félicien Akohoue <akohoue.f@gmail.com>
([ORCID](https://orcid.org/0000-0002-2160-0182))
